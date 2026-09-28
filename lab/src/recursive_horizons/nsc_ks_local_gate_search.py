"""Fail-closed trust-region campaign for the switched local incoming gate.

The numerical evaluator is injected through a narrow protocol.  This owner
controls admission, fresh-value acceptance, Jacobian checks and the immutable
event chain.  It does not turn a nodal search result into a gate certificate.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
import math
from pathlib import Path
import re
from typing import Protocol

import numpy as np

from .evidence_io import canonical_json_bytes, load_canonical_json
from .nsc_ks_evaluation_binding import write_bytes_atomic
from .nsc_ks_profile_identity import profile_identity
from .nsc_ks_retarded_trust_region import (
    RefreshPolicy,
    acceptance,
    broyden_update,
    flatten_constraints,
    levenberg_step,
    update_radius,
)
from .nsc_local_incoming_family import LocalIncomingFamily


COMPONENTS = (
    "field_space_time",
    "changed_history_UV_tail",
    "baseline_low_subgap",
    "upstream",
    "energy_interpolation",
    "covered_regions",
    "phase_value",
    "between_node",
    "arithmetic",
)
STATE_LAW = "C_Sigma[g]=U_g C_up U_g^dagger"
EVENT_RE = re.compile(r"^event-(\d{6})\.json$")
HEX64 = re.compile(r"^[0-9a-f]{64}$")


class SearchAdmissionError(ValueError):
    """The nine-component budget does not authorize a search campaign."""


class SearchBindingError(ValueError):
    """A checkpoint, evaluator receipt or dependency binding changed."""


class JacobianCheckError(ArithmeticError):
    """A full retarded Jacobian failed its centered-difference control."""


def _digest_bytes(raw: bytes) -> str:
    return sha256(raw).hexdigest()


def _digest_array(value: object) -> str:
    array = np.ascontiguousarray(np.asarray(value, dtype="<f8"))
    return _digest_bytes(array.tobytes())


def _finite_array(value: object, shape: tuple[int, ...], label: str) -> np.ndarray:
    array = np.array(value, dtype=float, copy=True)
    if array.shape != shape or not np.isfinite(array).all():
        raise ValueError(f"{label} must have shape {shape} and finite values")
    array.setflags(write=False)
    return array


def _hex64(value: str, label: str) -> str:
    if not isinstance(value, str) or HEX64.fullmatch(value) is None:
        raise ValueError(f"{label} must be a lowercase sha256 identity")
    return value


def _componentwise_upper_sum(bounds: list[np.ndarray]) -> np.ndarray:
    total = np.zeros(2, dtype=float)
    for bound in bounds:
        for index in range(2):
            total[index] = math.nextafter(total[index] + float(bound[index]), math.inf)
    return total


@dataclass(frozen=True, slots=True)
class Gate2Admission:
    budget_sha256: str
    component_sum: tuple[float, float]
    allocation_limit: float
    residual_reserve: float
    physical_tolerance: float

    def as_record(self) -> dict[str, object]:
        return {
            "budget_sha256": self.budget_sha256,
            "component_sum": list(self.component_sum),
            "allocation_limit": self.allocation_limit,
            "residual_reserve": self.residual_reserve,
            "physical_tolerance": self.physical_tolerance,
        }


def admit_gate2_budget(path: Path) -> Gate2Admission:
    """Require the complete v5 budget and at least 1e-11 residual headroom."""
    raw = Path(path).read_bytes()
    try:
        record = json.loads(raw)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise SearchAdmissionError("Gate-2 budget is not valid JSON") from error
    if record.get("schema") != "NSC-KS-GATE-ERROR-BUDGET-v5":
        raise SearchAdmissionError("Gate-2 requires NSC-KS-GATE-ERROR-BUDGET-v5")
    components = record.get("components")
    if not isinstance(components, dict) or set(components) != set(COMPONENTS):
        raise SearchAdmissionError("Gate-2 requires exactly the nine named components")
    if record.get("missing_components") not in ([], ()):
        raise SearchAdmissionError("Gate-2 budget still has missing components")
    if (record.get("covered_positive_families") != 60
            or record.get("covered_signed_families") != 120
            or record.get("all_declared_source_paths_closed") is not True):
        raise SearchAdmissionError("Gate-2 lacks complete 60/120 source coverage")
    bounds = []
    for name in COMPONENTS:
        row = components[name]
        if not isinstance(row, dict):
            raise SearchAdmissionError(f"{name} is not a component record")
        bound = np.asarray(row.get("bound"), dtype=float)
        if bound.shape != (2,) or not np.isfinite(bound).all() or np.any(bound < 0):
            raise SearchAdmissionError(f"{name} lacks a finite nonnegative N,beta bound")
        bounds.append(bound)
    total = _componentwise_upper_sum(bounds)
    declared = np.asarray(record.get("full_error_sum"), dtype=float)
    if (declared.shape != (2,) or not np.isfinite(declared).all()
            or np.any(declared < total)):
        raise SearchAdmissionError("full_error_sum does not enclose the nine components")
    limit = float(record.get("allocation_total", float("nan")))
    reserve = float(record.get("residual_reserve", float("nan")))
    tolerance = float(record.get("physical_tolerance", float("nan")))
    if (not np.isfinite([limit, reserve, tolerance]).all()
            or limit > 2e-11 or reserve < 1e-11 or tolerance != 3e-11
            or np.any(declared > 2e-11)
            or tolerance - float(np.max(declared)) < reserve):
        raise SearchAdmissionError("Gate-2 error and residual reserves do not fit 3e-11")
    if (record.get("constraint_order") != ["N", "beta"]
            or record.get("state_law") != STATE_LAW
            or record.get("numerical_indicators_used_as_bounds") is not False):
        raise SearchAdmissionError("Gate-2 scientific scope or bound provenance changed")
    status = record.get("status")
    if (not isinstance(status, str)
            or status.split(":", 1)[0] not in ("PASS", "CLOSED", "ENCLOSED")):
        raise SearchAdmissionError("Gate-2 budget status is not closed")
    return Gate2Admission(
        _digest_bytes(raw), tuple(map(float, declared)), limit, reserve, tolerance)


@dataclass(frozen=True, slots=True)
class SearchConfig:
    coefficient_count_per_function: int = 32
    node_count: int = 129
    nodal_target: float = 1e-11
    initial_trust_radius: float = 1e-2
    minimum_trust_radius: float = 1e-12
    maximum_trust_radius: float = 1.0
    finite_difference_relative_step: float = 2e-6
    finite_difference_absolute_step: float = 2e-8
    derivative_relative_tolerance: float = 2e-4
    derivative_absolute_tolerance: float = 2e-9
    max_trials: int = 50

    def __post_init__(self) -> None:
        if self.coefficient_count_per_function != 32 or self.node_count != 129:
            raise ValueError("production search is fixed to 32+32 coefficients and 129 nodes")
        finite_positive = (
            self.nodal_target,
            self.initial_trust_radius,
            self.minimum_trust_radius,
            self.maximum_trust_radius,
            self.finite_difference_relative_step,
            self.finite_difference_absolute_step,
            self.derivative_relative_tolerance,
            self.derivative_absolute_tolerance,
        )
        if not np.isfinite(finite_positive).all() or min(finite_positive) <= 0:
            raise ValueError("search tolerances and radii must be finite and positive")
        if not self.minimum_trust_radius <= self.initial_trust_radius <= self.maximum_trust_radius:
            raise ValueError("trust radii must be ordered")
        if (isinstance(self.max_trials, bool)
                or not isinstance(self.max_trials, (int, np.integer))
                or self.max_trials < 1):
            raise ValueError("max_trials must be a positive integer")

    @property
    def unknown_count(self) -> int:
        return 2 * self.coefficient_count_per_function

    def as_record(self) -> dict[str, object]:
        return {
            name: getattr(self, name)
            for name in self.__dataclass_fields__
        }


@dataclass(frozen=True, slots=True)
class SearchBinding:
    initial_profile_identity: str
    source_identity: str
    evaluator_identity: str
    evaluator_sha256: str
    solver_settings_sha256: str
    dependency_hashes: tuple[tuple[str, str], ...]
    state_law: str = STATE_LAW

    def __post_init__(self) -> None:
        _hex64(self.initial_profile_identity, "initial_profile_identity")
        _hex64(self.source_identity, "source_identity")
        _hex64(self.evaluator_identity, "evaluator_identity")
        _hex64(self.evaluator_sha256, "evaluator_sha256")
        _hex64(self.solver_settings_sha256, "solver_settings_sha256")
        if self.state_law != STATE_LAW:
            raise ValueError("switched state law binding changed")
        if (not self.dependency_hashes
                or len({name for name, _ in self.dependency_hashes}) != len(self.dependency_hashes)):
            raise ValueError("unique dependency hashes are required")
        for name, digest in self.dependency_hashes:
            if not isinstance(name, str) or not name or "\\" in name:
                raise ValueError("dependency names must be nonempty portable paths")
            _hex64(digest, "dependency digest")

    def as_record(self) -> dict[str, object]:
        return {
            "initial_profile_identity": self.initial_profile_identity,
            "source_identity": self.source_identity,
            "evaluator_identity": self.evaluator_identity,
            "evaluator_sha256": self.evaluator_sha256,
            "solver_settings_sha256": self.solver_settings_sha256,
            "dependency_hashes": {name: digest for name, digest in self.dependency_hashes},
            "state_law": self.state_law,
            "interval": "I=S(1)+[0.12,0.18]",
        }

    @property
    def digest(self) -> str:
        return _digest_bytes(canonical_json_bytes(self.as_record()))


@dataclass(frozen=True, slots=True)
class ValueEvaluation:
    gradient: object
    evaluation_identity: str
    coefficients_sha256: str
    evaluator_identity: str
    certified_uncertainty: float = 0.0
    fresh: bool = True

    def __post_init__(self) -> None:
        gradient = _finite_array(self.gradient, (129, 2), "value gradient")
        object.__setattr__(self, "gradient", gradient)
        for value, label in (
            (self.evaluation_identity, "evaluation_identity"),
            (self.coefficients_sha256, "coefficients_sha256"),
            (self.evaluator_identity, "evaluator_identity"),
        ):
            _hex64(value, label)
        uncertainty = float(self.certified_uncertainty)
        if not np.isfinite(uncertainty) or uncertainty < 0:
            raise ValueError("certified evaluation uncertainty must be finite and nonnegative")
        object.__setattr__(self, "certified_uncertainty", uncertainty)
        if self.fresh is not True:
            raise ValueError("search evaluations must be freshly computed")


@dataclass(frozen=True, slots=True)
class FullEvaluation(ValueEvaluation):
    jacobian: object = ()

    def __post_init__(self) -> None:
        super(FullEvaluation, self).__post_init__()
        jacobian = _finite_array(self.jacobian, (64, 129, 2), "full Jacobian")
        object.__setattr__(self, "jacobian", jacobian)


class SearchEvaluator(Protocol):
    evaluator_identity: str

    def preparation_supported(self, family: LocalIncomingFamily) -> bool: ...

    def value(self, family: LocalIncomingFamily, *, fresh: bool) -> ValueEvaluation: ...

    def full(self, family: LocalIncomingFamily, *, fresh: bool) -> FullEvaluation: ...


def coefficients_sha256(coefficients: object) -> str:
    array = _finite_array(coefficients, (2, 32), "history coefficients")
    return _digest_array(array)


def _require_evaluation(
    evaluation: ValueEvaluation,
    coefficients: np.ndarray,
    binding: SearchBinding,
) -> ValueEvaluation:
    if evaluation.coefficients_sha256 != coefficients_sha256(coefficients):
        raise SearchBindingError("evaluator receipt belongs to different coefficients")
    if evaluation.evaluator_identity != binding.evaluator_identity:
        raise SearchBindingError("evaluator identity changed")
    return evaluation


def selected_jacobian_columns(jacobian: object) -> tuple[int, ...]:
    matrix = np.asarray(jacobian, float)
    if matrix.shape != (64, 129, 2) or not np.isfinite(matrix).all():
        raise ValueError("full 64x129x2 Jacobian required")
    fixed = (0, 15, 31, 32, 47, 63)
    norms = np.linalg.norm(matrix.reshape(64, -1), axis=1)
    extra = [int(index) for index in np.argsort(norms)[::-1] if int(index) not in fixed][:2]
    return fixed + tuple(extra)


def _precheck(evaluator: SearchEvaluator, coefficients: np.ndarray) -> tuple[LocalIncomingFamily, bool, bool]:
    family = LocalIncomingFamily(coefficients)
    radius_positive = bool(family.radius_lower_bound() > 0)
    supported = bool(radius_positive and evaluator.preparation_supported(family))
    return family, radius_positive, supported


def verify_jacobian(
    evaluator: SearchEvaluator,
    coefficients: object,
    full: FullEvaluation,
    config: SearchConfig,
    binding: SearchBinding,
) -> dict[str, object]:
    """Check six fixed and two largest full-Jacobian columns by centered values."""
    coeff = _finite_array(coefficients, (2, 32), "history coefficients")
    _require_evaluation(full, coeff, binding)
    columns = selected_jacobian_columns(full.jacobian)
    rows = []
    for column in columns:
        base_value = float(coeff.ravel()[column])
        step = max(
            config.finite_difference_absolute_step,
            config.finite_difference_relative_step * max(1.0, abs(base_value)),
        )
        pair = None
        for _ in range(12):
            plus = np.array(coeff, copy=True)
            minus = np.array(coeff, copy=True)
            plus.ravel()[column] += step
            minus.ravel()[column] -= step
            plus_family, plus_radius, plus_support = _precheck(evaluator, plus)
            minus_family, minus_radius, minus_support = _precheck(evaluator, minus)
            if plus_radius and plus_support and minus_radius and minus_support:
                pair = (plus, minus, plus_family, minus_family)
                break
            step /= 2.0
        if pair is None:
            raise JacobianCheckError(f"column {column} has no admissible centered step")
        plus, minus, plus_family, minus_family = pair
        plus_eval = _require_evaluation(evaluator.value(plus_family, fresh=True), plus, binding)
        minus_eval = _require_evaluation(evaluator.value(minus_family, fresh=True), minus, binding)
        fd_gradient = (plus_eval.gradient - minus_eval.gradient) / (2.0 * step)
        _residual, fd_column = flatten_constraints(fd_gradient, fd_gradient[None, ...])
        exact_column = flatten_constraints(full.gradient, full.jacobian)[1][:, column]
        absolute_error = float(np.max(np.abs(fd_column[:, 0] - exact_column)))
        scale = float(max(
            np.max(np.abs(fd_column[:, 0])),
            np.max(np.abs(exact_column)),
            config.derivative_absolute_tolerance,
        ))
        uncertainty = (
            plus_eval.certified_uncertainty + minus_eval.certified_uncertainty
        ) / (2.0 * step)
        tolerance = (
            config.derivative_absolute_tolerance
            + config.derivative_relative_tolerance * scale
            + uncertainty
        )
        if absolute_error > tolerance:
            raise JacobianCheckError(
                f"Jacobian column {column} differs by {absolute_error} > {tolerance}")
        rows.append({
            "column": column,
            "centered_step": step,
            "absolute_error": absolute_error,
            "allowed_error": tolerance,
            "finite_difference_uncertainty": uncertainty,
            "plus_evaluation_identity": plus_eval.evaluation_identity,
            "minus_evaluation_identity": minus_eval.evaluation_identity,
        })
    receipt = {
        "full_evaluation_identity": full.evaluation_identity,
        "jacobian_sha256": _digest_array(full.jacobian),
        "selected_columns": list(columns),
        "checks": rows,
    }
    receipt["receipt_sha256"] = _digest_bytes(canonical_json_bytes(receipt))
    return receipt


@dataclass(frozen=True, slots=True)
class SearchOutcome:
    status: str
    coefficients: tuple[tuple[float, ...], tuple[float, ...]]
    component_maxima: tuple[float, float]
    merit: float
    trials: int
    accepted_trials: int
    last_event_sha256: str | None


def _manifest(
    initial: np.ndarray,
    metric: np.ndarray,
    config: SearchConfig,
    binding: SearchBinding,
    admission: Gate2Admission,
) -> dict[str, object]:
    if admission.budget_sha256 not in dict(binding.dependency_hashes).values():
        raise SearchBindingError("Gate-2 budget digest is absent from dependencies")
    return {
        "schema": "NSC-KS-LOCAL-GATE-SEARCH-MANIFEST-v1",
        "config": config.as_record(),
        "binding": binding.as_record(),
        "binding_sha256": binding.digest,
        "gate2_admission": admission.as_record(),
        "initial_coefficients": initial.tolist(),
        "initial_coefficients_sha256": coefficients_sha256(initial),
        "metric_sha256": _digest_array(metric),
        "prediction_is_not_a_residual": True,
        "physical_EXISTENCE_certificate": False,
        "physical_NONEXISTENCE_certificate": False,
    }


def _read_canonical(path: Path) -> dict[str, object]:
    raw = path.read_bytes()
    value = load_canonical_json(raw)
    if not isinstance(value, dict):
        raise SearchBindingError(f"{path.name} is not a JSON object")
    return value


def _open_journal(path: Path, expected: dict[str, object]) -> tuple[np.ndarray, float, int, int, str | None]:
    path.mkdir(parents=True, exist_ok=True)
    manifest_path = path / "manifest.json"
    raw = canonical_json_bytes(expected)
    if not manifest_path.exists():
        write_bytes_atomic(manifest_path, raw, exclusive=True)
    elif manifest_path.read_bytes() != raw:
        raise SearchBindingError("search manifest or dependency binding changed")
    events = sorted(
        (item for item in path.iterdir() if EVENT_RE.fullmatch(item.name)),
        key=lambda item: item.name,
    )
    coefficients = np.asarray(expected["initial_coefficients"], float)
    radius = float(expected["config"]["initial_trust_radius"])
    accepted = 0
    previous = None
    for expected_index, event_path in enumerate(events):
        match = EVENT_RE.fullmatch(event_path.name)
        if match is None or int(match.group(1)) != expected_index:
            raise SearchBindingError("event sequence is not contiguous")
        raw_event = event_path.read_bytes()
        event = load_canonical_json(raw_event)
        event_body = dict(event)
        declared_body_digest = event_body.pop("event_body_sha256", None)
        if declared_body_digest != _digest_bytes(canonical_json_bytes(event_body)):
            raise SearchBindingError("event body digest changed")
        if (event.get("schema") != "NSC-KS-LOCAL-GATE-SEARCH-EVENT-v1"
                or event.get("event_index") != expected_index
                or event.get("binding_sha256") != expected["binding_sha256"]
                or event.get("previous_event_sha256") != previous):
            raise SearchBindingError("event chain or binding changed")
        parent = np.asarray(event.get("parent_coefficients"), float)
        after = np.asarray(event.get("accepted_coefficients_after_event"), float)
        if (parent.shape != (2, 32) or after.shape != (2, 32)
                or not np.array_equal(parent, coefficients)):
            raise SearchBindingError("event does not continue the accepted state")
        if bool(event.get("accepted")):
            candidate = np.asarray(event.get("candidate_coefficients"), float)
            if candidate.shape != (2, 32) or not np.array_equal(candidate, after):
                raise SearchBindingError("accepted event coefficients changed")
            accepted += 1
        elif not np.array_equal(after, coefficients):
            raise SearchBindingError("rejected event changed accepted coefficients")
        coefficients = after
        radius = float(event.get("trust_radius_after"))
        if not np.isfinite(radius) or radius <= 0:
            raise SearchBindingError("event trust radius is invalid")
        previous = _digest_bytes(raw_event)
    return coefficients, radius, len(events), accepted, previous


def _write_event(journal: Path, event: dict[str, object]) -> str:
    index = int(event["event_index"])
    event = dict(event)
    event["event_body_sha256"] = _digest_bytes(canonical_json_bytes(event))
    raw = canonical_json_bytes(event)
    write_bytes_atomic(journal / f"event-{index:06d}.json", raw, exclusive=True)
    return _digest_bytes(raw)


def _component_maxima(gradient: np.ndarray) -> tuple[float, float]:
    return tuple(map(float, np.max(np.abs(gradient), axis=0)))


def run_search(
    evaluator: SearchEvaluator,
    initial_coefficients: object,
    metric: object,
    config: SearchConfig,
    binding: SearchBinding,
    admission: Gate2Admission,
    journal: Path,
    *,
    trial_limit: int | None = None,
) -> SearchOutcome:
    """Run or resume the admitted search, never certifying the physical gate."""
    initial = _finite_array(initial_coefficients, (2, 32), "initial coefficients")
    if profile_identity(LocalIncomingFamily(initial), include_normal_window=True) != binding.initial_profile_identity:
        raise SearchBindingError("initial history profile identity changed")
    metric_array = _finite_array(metric, (64, 64), "physical metric")
    eigenvalues = np.linalg.eigvalsh(metric_array)
    if eigenvalues[0] <= 0:
        raise ValueError("physical metric must be positive definite")
    if evaluator.evaluator_identity != binding.evaluator_identity:
        raise SearchBindingError("configured evaluator identity changed")
    expected = _manifest(initial, metric_array, config, binding, admission)
    coefficients, radius, event_index, accepted_count, previous = _open_journal(
        Path(journal), expected)
    family, radius_positive, supported = _precheck(evaluator, coefficients)
    if not radius_positive or not supported:
        raise SearchBindingError("resumed accepted state is outside radius/preparation support")

    full = _require_evaluation(evaluator.full(family, fresh=True), coefficients, binding)
    if not isinstance(full, FullEvaluation):
        raise TypeError("full evaluator must return FullEvaluation")
    jacobian_receipt = verify_jacobian(evaluator, coefficients, full, config, binding)
    gradient = np.asarray(full.gradient, float)
    _residual, jacobian = flatten_constraints(gradient, full.jacobian)
    current_uncertainty = full.certified_uncertainty
    policy = RefreshPolicy()
    trials_here = (
        config.max_trials - event_index
        if trial_limit is None else int(trial_limit)
    )
    if trials_here < 0 or event_index + trials_here > config.max_trials:
        raise ValueError("trial limit exceeds the configured campaign bound")

    maxima = _component_maxima(gradient)
    if max(maxima) <= config.nodal_target:
        return SearchOutcome(
            "NODAL_TARGET_REACHED; continuous certificate still required",
            tuple(tuple(map(float, row)) for row in coefficients), maxima,
            max(maxima), event_index, accepted_count, previous)

    completed = 0
    while completed < trials_here and event_index < config.max_trials:
        residual, _ = flatten_constraints(gradient, full.jacobian)
        proposal = levenberg_step(residual, jacobian, metric_array, radius)
        step = np.asarray(proposal["step"], float)
        candidate = np.asarray(coefficients, float).ravel() + step
        candidate = candidate.reshape(2, 32)
        candidate_family, candidate_radius, candidate_support = _precheck(evaluator, candidate)
        candidate_eval = None
        decision = None
        uncertainty_margin = None
        accepted = False
        if candidate_radius and candidate_support and proposal["predicted_merit"] < proposal["parent_merit"]:
            candidate_eval = _require_evaluation(
                evaluator.value(candidate_family, fresh=True), candidate, binding)
            decision = acceptance(
                gradient,
                candidate_eval.gradient,
                proposal["predicted_residual"],
                radius_positive=True,
                preparation_supported=True,
            )
            uncertainty_margin = current_uncertainty + candidate_eval.certified_uncertainty
            accepted = bool(
                decision["accepted"]
                and decision["measured_reduction"] > uncertainty_margin)
        ratio = None if decision is None or not np.isfinite(decision["prediction_ratio"]) else float(
            decision["prediction_ratio"])
        policy_before = policy
        policy, refresh = policy.update(accepted=accepted, prediction_ratio=(
            float("-inf") if ratio is None else ratio))
        radius_after = update_radius(
            radius,
            float("-inf") if ratio is None else ratio,
            accepted,
            minimum=config.minimum_trust_radius,
            maximum=config.maximum_trust_radius,
        )
        after = np.array(candidate if accepted else coefficients, copy=True)
        event = {
            "schema": "NSC-KS-LOCAL-GATE-SEARCH-EVENT-v1",
            "event_index": event_index,
            "binding_sha256": binding.digest,
            "previous_event_sha256": previous,
            "accepted": accepted,
            "parent_coefficients": coefficients.tolist(),
            "candidate_coefficients": candidate.tolist(),
            "accepted_coefficients_after_event": after.tolist(),
            "trust_radius_before": radius,
            "trust_radius_after": radius_after,
            "step_metric_norm": float(proposal["metric_norm"]),
            "parent_merit": float(proposal["parent_merit"]),
            "predicted_merit": float(proposal["predicted_merit"]),
            "measured_merit": None if candidate_eval is None else float(
                np.max(np.abs(candidate_eval.gradient))),
            "predicted_reduction": None if decision is None else float(
                decision["predicted_reduction"]),
            "measured_reduction": None if decision is None else float(
                decision["measured_reduction"]),
            "prediction_ratio": ratio,
            "confirmed_reduction_margin": uncertainty_margin,
            "radius_positive": candidate_radius,
            "preparation_supported": candidate_support,
            "fresh_evaluation_identity": None if candidate_eval is None else candidate_eval.evaluation_identity,
            "fresh_evaluation_coefficients_sha256": None if candidate_eval is None else candidate_eval.coefficients_sha256,
            "jacobian_evaluation_identity": full.evaluation_identity,
            "jacobian_sha256": _digest_array(jacobian),
            "jacobian_check": jacobian_receipt,
            "broyden_counters_before": {
                "accepted_since_refresh": policy_before.accepted_since_refresh,
                "poor_ratios_since_refresh": policy_before.poor_ratios_since_refresh,
            },
            "broyden_counters_after": {
                "accepted_since_refresh": policy.accepted_since_refresh,
                "poor_ratios_since_refresh": policy.poor_ratios_since_refresh,
            },
            "refresh_required_after_event": refresh,
            "prediction_is_not_a_residual": True,
            "physical_EXISTENCE_certificate": False,
            "physical_NONEXISTENCE_certificate": False,
        }
        previous = _write_event(Path(journal), event)
        event_index += 1
        completed += 1
        radius = radius_after
        if accepted:
            change = np.concatenate((
                candidate_eval.gradient[:, 0] - gradient[:, 0],
                candidate_eval.gradient[:, 1] - gradient[:, 1],
            ))
            jacobian = broyden_update(jacobian, step, change)
            coefficients = after
            gradient = np.asarray(candidate_eval.gradient, float)
            current_uncertainty = candidate_eval.certified_uncertainty
            accepted_count += 1
            maxima = _component_maxima(gradient)
            if max(maxima) <= config.nodal_target:
                break
        if refresh:
            family, radius_positive, supported = _precheck(evaluator, coefficients)
            if not radius_positive or not supported:
                raise SearchBindingError("accepted state failed pre-refresh support")
            full = _require_evaluation(evaluator.full(family, fresh=True), coefficients, binding)
            difference = float(np.max(np.abs(full.gradient - gradient)))
            if difference > full.certified_uncertainty + current_uncertainty:
                raise SearchBindingError("fresh full refresh disagrees with the accepted value")
            gradient = np.asarray(full.gradient, float)
            current_uncertainty = full.certified_uncertainty
            _residual, jacobian = flatten_constraints(gradient, full.jacobian)
            jacobian_receipt = verify_jacobian(evaluator, coefficients, full, config, binding)

    maxima = _component_maxima(gradient)
    status = (
        "NODAL_TARGET_REACHED; continuous certificate still required"
        if max(maxima) <= config.nodal_target
        else "OPEN: admitted search incomplete or nodal target not reached"
    )
    return SearchOutcome(
        status,
        tuple(tuple(map(float, row)) for row in coefficients),
        maxima,
        max(maxima),
        event_index,
        accepted_count,
        previous,
    )
