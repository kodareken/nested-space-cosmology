#!/usr/bin/env python3
"""Run the admitted 129-node, 64-direction switched-state gate search.

The command refuses to build a scientific evaluator until the successor v5
budget closes all nine components and the full 60/120 source inventory.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import sys

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import derive_nsc_ks_coupled_retained_control as Retained
import derive_nsc_ks_gate_value as Value
from recursive_horizons.evidence_io import canonical_json_bytes
from recursive_horizons.nsc_dirac_source_phase import (
    _normal_coordinate, formal_source_phase_coefficient)
from recursive_horizons.nsc_evolved_incoming_constraints import (
    compatible_history_slots, surface_geometry_response)
from recursive_horizons.nsc_ks_cutoff_bridge_evaluator import (
    signed_family_angular_square_weights, source_edge_from_phase)
from recursive_horizons.nsc_ks_history_evaluator import KSHistoryEvaluator
from recursive_horizons.nsc_ks_local_gate_search import (
    FullEvaluation,
    SearchAdmissionError,
    SearchBinding,
    SearchConfig,
    ValueEvaluation,
    admit_gate2_budget,
    coefficients_sha256,
    run_search,
)
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_ks_retarded_trust_region import physical_metric
from recursive_horizons.nsc_ks_source_envelope import computational_z_grid
from recursive_horizons.nsc_ks_evaluation_binding import write_bytes_atomic
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily


BUDGET = ROOT / "results/development/nsc-ks-gate-budget-v5.json"
HISTORY = ROOT / "results/development/nsc-ks-gate-history-lm-broyden.json"
PRINCIPAL = ROOT / "results/development/nsc-incoming-surface-principal.json"
ITERATE5 = ROOT / "results/development/nsc-ks-coupled-newton-n32-truncated-iterate5.json"
ITERATE6 = ROOT / "results/development/nsc-ks-coupled-newton-n32-truncated-iterate6.json"
JOURNAL = ROOT / "results/development/artifacts/nsc-ks-local-gate-search-v1"
SOLVER = {
    "grid_nodes": 1024,
    "length": 0.4,
    "family_energy_intervals": [[0.0, 160.0], [0.0, 320.0]],
    "interpolant_degree": 48,
    "rtol": 5e-14,
    "atol": 5e-19,
    "max_step": 1 / 8192,
    "integrator": "dop853",
    "step_control": "joint",
    "phase_nodes": 192,
    "target_nodes": 129,
    "retarded_directions": 64,
}
VALUE_WORKERS = 9
VALUE_CPU_BUDGET = 14400.0


def digest(path: Path) -> str:
    return sha256(Path(path).read_bytes()).hexdigest()


def _receipt_identity(kind, evaluator_identity, coefficients, gradient, jacobian=None):
    record = {
        "kind": kind,
        "evaluator_identity": evaluator_identity,
        "coefficients_sha256": coefficients_sha256(coefficients),
        "gradient_sha256": sha256(
            np.ascontiguousarray(gradient, dtype="<f8").tobytes()).hexdigest(),
        "jacobian_sha256": None if jacobian is None else sha256(
            np.ascontiguousarray(jacobian, dtype="<f8").tobytes()).hexdigest(),
    }
    return sha256(canonical_json_bytes(record)).hexdigest()


class ProductionEvaluator:
    """Rebuild the complete evaluator for every receipt, preventing cache reuse."""

    def __init__(self, context, admission, evaluator_identity):
        self._context = context
        self._admission = admission
        self.evaluator_identity = evaluator_identity
        self._entries = tuple(
            entry
            for key in sorted(context["archived"])
            for entry in context["archive"].family_entries(key)
        )
        if len(context["archived"]) != 60:
            raise ValueError("production evaluator requires all 60 positive families")
        self._rho_up = {float(batch.rho_up) for batch, _channel in self._entries}
        if len(self._rho_up) != 1:
            raise ValueError("fixed upstream source must have one rho_up")
        self._nodes = LocalIncomingFamily(
            np.asarray(json.loads(HISTORY.read_text())["history"]["coefficients"], float)
        ).collocation_nodes(129)
        self._grid = computational_z_grid(SOLVER["grid_nodes"], SOLVER["length"])

    def preparation_supported(self, family):
        rho_up = next(iter(self._rho_up))
        return bool(
            family.normal_outer == 0.03
            and family.normal_inner == 0.007
            and _normal_coordinate(rho_up) <= -family.normal_outer
        )

    def _value_history(self, family):
        coefficients = np.asarray(family.coefficients, float)
        coefficient_digest = coefficients_sha256(coefficients)
        label = "search-v1-" + coefficient_digest[:20]
        history_path = ROOT / "results/development" / f"nsc-ks-gate-history-{label}.json"
        record = {
            "schema": "NSC-KS-GATE-HISTORY-v1",
            "profile_identity": profile_identity(family, include_normal_window=True),
            "history": family.description(),
            "prediction_is_not_a_residual": True,
        }
        raw = json.dumps(record, sort_keys=True, indent=2, allow_nan=False).encode() + b"\n"
        if history_path.exists():
            if history_path.read_bytes() != raw:
                raise ValueError("value-only history binding changed")
        else:
            write_bytes_atomic(history_path, raw, exclusive=True)
        output = Value.output_paths(label)[1]
        settings = {
            key: SOLVER[key]
            for key in ("grid_nodes", "interpolant_degree", "rtol", "atol",
                        "max_step", "phase_nodes", "length", "integrator", "step_control")
        }
        settings["degree"] = settings.pop("interpolant_degree")
        if not output.exists():
            result = Value.run(
                129, VALUE_WORKERS, VALUE_CPU_BUDGET, label, settings,
                str(history_path))
            if result["completed_positive_families"] != result["required_positive_families"]:
                raise RuntimeError("value-only evaluation exhausted its CPU budget; resume the same journal")
        gradient = Value.reconstruct_gradient(label)
        top = json.loads(output.read_text())
        return gradient, {
            "record_sha256": digest(output),
            "evaluation_identity": top["evaluation_identity"],
        }

    def _full_evaluate(self, family):
        matter = np.zeros((129, 2), float)
        tangent = np.zeros((64, 129, 2), float)
        family_records = {}
        for key in sorted(self._context["archived"]):
            entries = self._context["archive"].family_entries(key)
            raw = KSHistoryEvaluator(
                entries,
                baseline_gradient=self._context["baseline"],
                coefficients=self._context["coeff"],
                energy_interval=Retained.interpolation_interval(
                    self._context["archived"][key]),
                interpolant_degree=SOLVER["interpolant_degree"],
                z_grid=self._grid,
                solve_nodes=self._nodes,
                verification_nodes=self._nodes,
                rtol=SOLVER["rtol"],
                atol=SOLVER["atol"],
                max_step=SOLVER["max_step"],
                integrator=SOLVER["integrator"],
                step_control=SOLVER["step_control"],
            ).evaluate(family, self._nodes)
            matter += sum(np.asarray(value, float)
                          for value in raw["family_corrections"].values())
            tangent += sum(np.asarray(value, float)
                           for value in raw["family_matter_tangents"].values())
            overlap = set(family_records) & set(raw["family_records"])
            if overlap:
                raise ValueError("signed source family counted twice")
            family_records.update(raw["family_records"])
        if len(family_records) != 120:
            raise ValueError("full evaluator did not retain all 120 signed contributions")
        slots, slot_tangents = compatible_history_slots(family.metric(), self._nodes, 64)
        geometry = surface_geometry_response(
            slots, slot_tangents, self._context["coeff"])
        baseline = np.broadcast_to(self._context["baseline"], (129, 2)).copy()
        rho_up = next(iter(self._rho_up))
        phase = formal_source_phase_coefficient(
            family, self._nodes, angular=1.0, rho_up=rho_up,
            gauss_nodes=SOLVER["phase_nodes"])
        weights = signed_family_angular_square_weights(family_records)
        edge, edge_tangent = source_edge_from_phase(
            phase, weights, self._context["coeff"]["a"])
        return (
            baseline + geometry["action_gradient_change"] + matter + edge,
            geometry["action_gradient_tangent"] + tangent + edge_tangent,
        )

    def value(self, family, *, fresh):
        if fresh is not True:
            raise ValueError("production value evaluation must be fresh")
        coefficients = np.asarray(family.coefficients, float)
        gradient, evidence = self._value_history(family)
        identity = sha256(canonical_json_bytes({
            "receipt": _receipt_identity(
                "value", self.evaluator_identity, coefficients, gradient),
            "evidence": evidence,
        })).hexdigest()
        return ValueEvaluation(
            gradient,
            identity,
            coefficients_sha256(coefficients),
            self.evaluator_identity,
            max(self._admission.component_sum),
            True,
        )

    def full(self, family, *, fresh):
        if fresh is not True:
            raise ValueError("production full evaluation must be fresh")
        coefficients = np.asarray(family.coefficients, float)
        gradient, jacobian = self._full_evaluate(family)
        identity = _receipt_identity(
            "full", self.evaluator_identity, coefficients, gradient, jacobian)
        return FullEvaluation(
            gradient,
            identity,
            coefficients_sha256(coefficients),
            self.evaluator_identity,
            max(self._admission.component_sum),
            True,
            jacobian,
        )


def initial_family():
    record = json.loads(HISTORY.read_text())
    family = LocalIncomingFamily(np.asarray(record["history"]["coefficients"], float))
    if profile_identity(family, include_normal_window=True) != record["profile_identity"]:
        raise ValueError("initial history identity changed")
    return family


def principal_matrix():
    intervals = json.loads(PRINCIPAL.read_text())["matrix"]["matrix_intervals"]
    matrix = np.asarray([[0.5 * (lo + hi) for lo, hi in row] for row in intervals], float)
    if matrix.shape != (2, 2) or not np.isfinite(matrix).all():
        raise ValueError("certified principal midpoint must be a finite 2x2 matrix")
    return matrix


def initial_radius(metric):
    parent = np.asarray(json.loads(ITERATE5.read_text())["history"]["coefficients"], float)
    child = np.asarray(json.loads(ITERATE6.read_text())["history"]["coefficients"], float)
    step = (child - parent).ravel()
    radius = float(np.sqrt(step @ metric @ step))
    if not np.isfinite(radius) or radius <= 0:
        raise ValueError("recorded accepted repair does not define a trust radius")
    return radius


def build_campaign(admission, budget_path=BUDGET):
    budget_path = Path(budget_path).resolve()
    try:
        budget_path.relative_to(ROOT.resolve())
    except ValueError as error:
        raise ValueError("production budget must be inside the repository") from error
    context = Retained.context()
    family = initial_family()
    metric = physical_metric(family, principal_matrix())
    radius = initial_radius(metric)
    solver_hash = sha256(canonical_json_bytes(SOLVER)).hexdigest()
    source_record = {
        "archive_input_hashes": dict(sorted(context["archive"].input_hashes.items())),
        "positive_families": [list(key) for key in sorted(context["archived"])],
        "analytic_zero_groups": [int(row["index"]) for row in context["ell0"]],
    }
    source_identity = sha256(canonical_json_bytes(source_record)).hexdigest()
    evaluator_record = {
        "source_identity": source_identity,
        "solver": SOLVER,
        "state_law": "C_Sigma[g]=U_g C_up U_g^dagger",
    }
    evaluator_identity = sha256(canonical_json_bytes(evaluator_record)).hexdigest()
    owned = (
        Path(__file__),
        ROOT / "src/recursive_horizons/nsc_ks_local_gate_search.py",
        ROOT / "src/recursive_horizons/nsc_ks_retarded_trust_region.py",
        ROOT / "src/recursive_horizons/nsc_ks_history_evaluator.py",
        ROOT / "src/recursive_horizons/nsc_ks_cutoff_bridge_evaluator.py",
        HISTORY,
        PRINCIPAL,
        ITERATE5,
        ITERATE6,
        budget_path,
    )
    hashes = {
        str(path.relative_to(ROOT)).replace("\\", "/"): digest(path)
        for path in owned
    }
    hashes.update({str(name).replace("\\", "/"): value
                   for name, value in context["archive"].input_hashes.items()})
    binding = SearchBinding(
        profile_identity(family, include_normal_window=True),
        source_identity,
        evaluator_identity,
        digest(Path(__file__)),
        solver_hash,
        tuple(sorted(hashes.items())),
    )
    config = SearchConfig(
        initial_trust_radius=radius,
        maximum_trust_radius=8 * radius,
    )
    evaluator = ProductionEvaluator(context, admission, evaluator_identity)
    return evaluator, np.asarray(family.coefficients, float), metric, config, binding


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--audit-admission", action="store_true")
    mode.add_argument("--run", action="store_true")
    parser.add_argument("--budget", type=Path, default=BUDGET)
    parser.add_argument("--journal", type=Path, default=JOURNAL)
    parser.add_argument("--trial-limit", type=int)
    args = parser.parse_args(argv)
    try:
        admission = admit_gate2_budget(args.budget)
    except (FileNotFoundError, SearchAdmissionError) as error:
        print(json.dumps({
            "status": "OPEN: Gate-2 admission refused",
            "budget": str(args.budget),
            "reason": str(error),
            "search_started": False,
        }, indent=2))
        return 2
    if args.audit_admission:
        print(json.dumps({
            "status": "PASS: production search admitted",
            "gate2": admission.as_record(),
            "search_started": False,
        }, indent=2))
        return 0
    evaluator, initial, metric, config, binding = build_campaign(admission, args.budget)
    outcome = run_search(
        evaluator, initial, metric, config, binding, admission, args.journal,
        trial_limit=args.trial_limit)
    print(json.dumps({
        "status": outcome.status,
        "component_maxima": list(outcome.component_maxima),
        "merit": outcome.merit,
        "trials": outcome.trials,
        "accepted_trials": outcome.accepted_trials,
        "last_event_sha256": outcome.last_event_sha256,
        "physical_EXISTENCE_certificate": False,
        "physical_NONEXISTENCE_certificate": False,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
