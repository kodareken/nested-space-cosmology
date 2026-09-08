#!/usr/bin/env python3
"""Recompute the curated result graph without modifying the checkout."""

from __future__ import annotations

import argparse
import ast
from concurrent.futures import FIRST_COMPLETED, Future, ThreadPoolExecutor, wait
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "results" / "manifest.json"
CORE_FILES = (
    "src/recursive_horizons/__init__.py",
    "src/recursive_horizons/evidence_io.py",
    "src/recursive_horizons/unified_action.py",
    "src/recursive_horizons/nsc_regulated.py",
    "src/recursive_horizons/nsc_geometric_chain.py",
)


class ReproductionError(RuntimeError):
    """The public result graph could not be reproduced."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("portable", "exact"), default="portable")
    parser.add_argument("--jobs", default="auto", help="positive integer or 'auto'")
    parser.add_argument("--keep-workspace", action="store_true")
    parser.add_argument("--json-summary", type=Path)
    return parser.parse_args()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_json(value: object) -> bytes:
    return json.dumps(
        value,
        allow_nan=False,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("ascii")


def resolve_jobs(value: str) -> int:
    if value == "auto":
        return min(8, max(1, os.cpu_count() or 1))
    try:
        jobs = int(value)
    except ValueError as exc:
        raise ReproductionError("--jobs must be a positive integer or 'auto'") from exc
    if jobs < 1:
        raise ReproductionError("--jobs must be positive")
    return jobs


def load_manifest() -> dict[str, Any]:
    value = json.loads(MANIFEST.read_text(encoding="utf-8"))
    if value.get("schema") != "NSC-PUBLIC-RESULT-MANIFEST-v1":
        raise ReproductionError("unexpected result manifest schema")
    steps = value.get("steps")
    historical = value.get("historical_result_count", 58)
    spec_path = ROOT / value.get("release_spec", "results/release-spec.json")
    spec = json.loads(spec_path.read_text(encoding="utf-8"))
    expected_count = 58 + len(spec["scoped_follow_ups"])
    if (not isinstance(steps, list) or len(steps) != expected_count
            or value.get("result_count") != expected_count):
        raise ReproductionError("public result count differs from the explicit release specification")
    if sha256(spec_path) != value.get("release_spec_sha256"):
        raise ReproductionError("release specification hash mismatch")
    if historical != 58:
        raise ReproductionError("historical result count must remain 58")
    return value


def validate_checkout(manifest: dict[str, Any]) -> None:
    spec_path = ROOT / manifest["release_spec"]
    if sha256(spec_path) != manifest["release_spec_sha256"]:
        raise ReproductionError("release specification hash mismatch")
    spec = json.loads(spec_path.read_text(encoding="utf-8"))
    scoped_specs = {row["output"]: row for row in spec["scoped_follow_ups"]}
    if (manifest.get("historical_result_count") != 58
            or len(manifest["steps"]) != 58 + len(scoped_specs)
            or {row["output"] for row in manifest["steps"][58:]} != set(scoped_specs)):
        raise ReproductionError("manifest graph differs from the historical and scoped release specification")
    for step in manifest["steps"][58:]:
        declared = scoped_specs[step["output"]]
        for key in ("artifact_id", "category"):
            if key in declared and step.get(key) != declared[key]:
                raise ReproductionError(f"manifest {key} differs from release specification: {step['output']}")
        for key in ("generator", "generator_args", "json_format", "comparison_policy",
                    "identity_policy", "dependencies", "source_dependencies", "auxiliary_inputs",
                    "follow_up_source_commit"):
            if step.get(key) != declared.get(key):
                raise ReproductionError(f"manifest {key} differs from release specification: {step['output']}")
    for item in (spec["import_files"] + spec["retained_byte_identical_files"]
                 + spec["preserved_64_scientific_files"]
                 + spec.get("preserved_75_scientific_files", [])
                 + spec.get("preserved_77_scientific_files", []) + [spec["retained_comparator"]]):
        relative = item["path"]
        if not (ROOT / relative).is_file() or sha256(ROOT / relative) != item["sha256"]:
            raise ReproductionError(f"pinned release input mismatch: {relative}")
    outputs: set[str] = set()
    for step in manifest["steps"]:
        output = str(step["output"])
        generator = str(step["generator"])
        if output in outputs:
            raise ReproductionError(f"duplicate result output: {output}")
        if output in set(step["dependencies"]):
            raise ReproductionError(f"self dependency: {output}")
        for relative, expected in (
            (output, step["output_sha256"]),
            (generator, step["generator_sha256"]),
        ):
            path = ROOT / relative
            if not path.is_file():
                raise ReproductionError(f"declared public file is absent: {relative}")
            observed = sha256(path)
            if observed != expected:
                raise ReproductionError(
                    f"tracked hash mismatch for {relative}: {observed} != {expected}"
                )
        if not set(step["dependencies"]) <= outputs:
            missing = sorted(set(step["dependencies"]) - outputs)
            raise ReproductionError(
                f"manifest is not topological at {output}: {', '.join(missing)}"
            )
        outputs.add(output)
        value = json.loads((ROOT / output).read_text(encoding="utf-8"))
        validate_identity(value, step)
        validate_authenticated_inputs(ROOT, value, step, use_auxiliary=False)
        if set(step.get("source_dependency_hashes", {})) != set(step.get("source_dependencies", [])) and "source_dependency_hashes" in step:
            raise ReproductionError(f"incomplete source hash closure: {output}")
        for relative, expected in step.get("source_dependency_hashes", {}).items():
            if not (ROOT / relative).is_file() or sha256(ROOT / relative) != expected:
                raise ReproductionError(f"source dependency mismatch: {relative}")
        for auxiliary in step.get("auxiliary_inputs", []):
            relative = str(auxiliary["path"])
            expected = str(auxiliary["sha256"])
            path = ROOT / relative
            if not path.is_file() or sha256(path) != expected:
                raise ReproductionError(
                    f"auxiliary provenance input mismatch: {relative}"
                )
    if manifest["frontier_output"] not in outputs:
        raise ReproductionError("frontier output is not present in the graph")
    follow_up = manifest.get("follow_up_output")
    if follow_up is not None and follow_up not in outputs:
        raise ReproductionError("follow-up output is not present in the graph")


def copy_relative(source_root: Path, target_root: Path, relative: str) -> None:
    source = source_root / relative
    if not source.is_file():
        raise ReproductionError(f"copy source is absent: {relative}")
    target = target_root / relative
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)


def prepare_workspace(path: Path, manifest: dict[str, Any]) -> Path:
    work = path / "work"
    expected = path / "expected"
    for relative in CORE_FILES:
        copy_relative(ROOT, work, relative)
    scripts: set[str] = set()
    for step in manifest["steps"]:
        scripts.add(str(step["generator"]))
        scripts.update(str(item) for item in step.get("source_dependencies", []))
        copy_relative(ROOT, expected, str(step["output"]))
        auxiliary_root = Path("provenance") / Path(str(step["output"])).stem
        for auxiliary in step.get("auxiliary_inputs", []):
            relative = str(auxiliary["path"])
            source = ROOT / relative
            if not source.is_file():
                raise ReproductionError(f"auxiliary provenance input is absent: {relative}")
            target = work / auxiliary_root / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
    for relative in sorted(scripts):
        copy_relative(ROOT, work, relative)
    (work / "results").mkdir(parents=True, exist_ok=True)
    return work


def run_generator(work: Path, step: dict[str, Any]) -> tuple[dict[str, Any], float]:
    environment = os.environ.copy()
    environment.update(
        {
            "PYTHONHASHSEED": "0",
            "OPENBLAS_NUM_THREADS": "1",
            "OMP_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "VECLIB_MAXIMUM_THREADS": "1",
        }
    )
    started = time.monotonic()
    command = [sys.executable, str(work / step["generator"])]
    command.extend(str(item) for item in step.get("generator_args", []))
    process = subprocess.run(
        command,
        cwd=work,
        env=environment,
        capture_output=True,
        text=True,
    )
    elapsed = time.monotonic() - started
    if process.returncode:
        detail = (process.stderr or process.stdout)[-4000:]
        raise ReproductionError(
            f"generator failed for {step['artifact_id']} ({process.returncode}):\n{detail}"
        )
    output = work / step["output"]
    if not output.is_file():
        raise ReproductionError(f"generator published no output: {step['output']}")
    raw = output.read_bytes()
    value = json.loads(raw)
    if step.get("json_format") != "pretty" and raw != canonical_json(value):
        raise ReproductionError(f"result is not canonical JSON: {step['output']}")
    validate_identity(value, step)
    return value, elapsed


def validate_identity(value: dict[str, Any], step: dict[str, Any]) -> None:
    policy = step.get("identity_policy", {"kind": "terminal"})
    if "schema" in policy and value.get("schema") != policy["schema"]:
        raise ReproductionError(f"invalid result schema: {step['output']}")
    if policy["kind"] == "schema_gate":
        if "artifact_id" in value or "terminal" in value:
            raise ReproductionError(f"unexpected identity fields in schema-gated result: {step['output']}")
        if value.get("classification") != policy["classification"]:
            raise ReproductionError(f"invalid schema-gated classification: {step['output']}")
        gate = value.get("gate")
        if (not isinstance(gate, dict) or set(gate) != set(policy["gate"])
                or any(gate[key] is not expected for key, expected in policy["gate"].items())):
            raise ReproductionError(f"invalid schema-gated completion: {step['output']}")
        return
    if value.get("artifact_id") != step["artifact_id"]:
        raise ReproductionError(f"invalid result identity: {step['output']}")
    if policy["kind"] == "plateau_scope":
        if (value.get("scope") != policy["scope"]
                or value.get("observational_audit", {}).get("prediction_claim") is not False):
            raise ReproductionError(f"invalid plateau scope or prediction claim: {step['output']}")
    elif policy["kind"] != "terminal" or value.get("terminal") is not True:
        raise ReproductionError(f"invalid terminal policy: {step['output']}")


def hash_bindings(value: dict[str, Any]):
    """Yield exact JSON hash locations and file references for supported schemas."""
    for key, collection in value.items():
        if key == "source_hashes" and isinstance(collection, list):
            for index, entry in enumerate(collection):
                if not isinstance(entry, dict) or not {"path", "sha256"} <= set(entry):
                    raise ReproductionError("malformed source_hashes list entry")
                yield (key, index, "sha256"), entry["path"], entry["sha256"], entry.get("artifact_id")
        elif key in {"source_hashes", "input_hashes", "authenticated_input_hashes"}:
            if not isinstance(collection, dict):
                raise ReproductionError(f"malformed {key} mapping")
            for relative, expected in collection.items():
                if not isinstance(relative, str) or not isinstance(expected, str):
                    raise ReproductionError(f"malformed {key} hash entry")
                yield (key, relative), relative, expected, None
        elif key.startswith("authenticated"):
            if isinstance(collection, dict) and {"path", "sha256"} <= set(collection):
                entries = [((key,), collection)]
            elif isinstance(collection, dict):
                entries = [((key, item), entry) for item, entry in collection.items()]
            elif isinstance(collection, list):
                entries = [((key, index), entry) for index, entry in enumerate(collection)]
            else:
                raise ReproductionError(f"malformed {key} collection")
            for location, entry in entries:
                if not isinstance(entry, dict) or not {"path", "sha256"} <= set(entry):
                    raise ReproductionError(f"malformed entry in {key}")
                yield location + ("sha256",), entry["path"], entry["sha256"], entry.get("artifact_id")
        elif key == "independent_tetrad_record_sha256":
            yield (key,), "results/nsc-4-dirac-tetrad.json", collection, "NSC-4-DIRAC-TETRAD"


def authenticated_entries(value: dict[str, Any]):
    """Compatibility view; mutation uses hash_bindings' exact JSON locations."""
    for location, relative, expected, artifact in hash_bindings(value):
        if str(location[0]).startswith("authenticated"):
            yield {"path": relative, "sha256": expected, "artifact_id": artifact}


def auxiliary_lookup(work: Path, step: dict[str, Any] | None) -> dict[str, Path]:
    if not step:
        return {}
    root = work / "provenance" / Path(str(step["output"])).stem
    return {str(item["path"]): root / str(item["path"])
            for item in step.get("auxiliary_inputs", [])}


def validate_authenticated_inputs(
    work: Path, value: dict[str, Any], step: dict[str, Any] | None = None,
    *, use_auxiliary: bool = True,
) -> set[tuple]:
    lookup = auxiliary_lookup(work, step) if use_auxiliary else {}
    auxiliary_paths = {str(item["path"]) for item in (step or {}).get("auxiliary_inputs", [])}
    generated = set((step or {}).get("dependencies", []))
    validated_dynamic: set[tuple] = set()
    for location, relative, expected, artifact_id in hash_bindings(value):
        if (not isinstance(relative, str) or not isinstance(expected, str)
                or len(expected) != 64 or any(char not in "0123456789abcdef" for char in expected)
                or Path(relative).is_absolute() or ".." in Path(relative).parts):
            raise ReproductionError("malformed authenticated file or hash")
        if (step is not None and relative.startswith("results/")
                and relative not in generated | auxiliary_paths):
            raise ReproductionError(f"authenticated result is missing from dependency graph: {relative}")
        path = lookup.get(relative, work / relative)
        if not path.is_file() or sha256(path) != expected:
            raise ReproductionError(f"authenticated input mismatch: {relative}")
        if artifact_id is not None and relative.startswith("results/"):
            actual_id = json.loads(path.read_text(encoding="utf-8")).get("artifact_id")
            if actual_id != artifact_id:
                raise ReproductionError(f"authenticated artifact mismatch: {relative}")
        if (relative in generated and relative not in auxiliary_paths
                and relative.startswith("results/") and location[0] != "source_hashes"):
            validated_dynamic.add(location)
    if "source_sha256" in value and step is not None:
        if sha256(work / step["generator"]) != value["source_sha256"]:
            raise ReproductionError("generator source_sha256 mismatch")
    if "comparison_function_sha256" in value:
        source = (work / "scripts/check_nsc_scale_closure.py").read_text(encoding="utf-8")
        node = next(node for node in ast.parse(source).body
                    if isinstance(node, ast.FunctionDef) and node.name == "compare")
        body = "".join(source.splitlines(keepends=True)[node.lineno-1:node.end_lineno])
        if hashlib.sha256(body.encode()).hexdigest() != value["comparison_function_sha256"]:
            raise ReproductionError("comparison function hash mismatch")
    return validated_dynamic


def normalize_dynamic_hashes(value: dict[str, Any], validated_locations: set[tuple]) -> dict[str, Any]:
    """Normalize only generated-result bindings already authenticated by the caller."""
    normalized = json.loads(json.dumps(value))
    allowed = {location for location, relative, _, _ in hash_bindings(value)
               if relative.startswith("results/") and location[0] != "source_hashes"}
    if not validated_locations <= allowed:
        raise ReproductionError("normalization requested for a non-result hash")
    for location in validated_locations:
        parent = normalized
        for part in location[:-1]:
            parent = parent[part]
        parent[location[-1]] = "<validated-generated-result>"
    return normalized


def compare_portable(
    expected: Any,
    actual: Any,
    *,
    path: str,
    relative_tolerance: float,
    absolute_tolerance: float,
    compare_numbers: bool,
    numeric_overrides: dict[str, tuple[float, float]] | None = None,
) -> None:
    if numeric_overrides and path in numeric_overrides:
        relative_tolerance, absolute_tolerance = numeric_overrides[path]
    if isinstance(expected, bool) or isinstance(actual, bool):
        if expected is not actual:
            raise ReproductionError(f"portable mismatch at {path}: {actual!r} != {expected!r}")
        return
    if isinstance(expected, (int, float)) and isinstance(actual, (int, float)):
        if not compare_numbers:
            return
        if isinstance(expected, int):
            if not isinstance(actual, int) or expected != actual:
                raise ReproductionError(f"integer mismatch at {path}: {actual!r} != {expected!r}")
            return
        if not math.isclose(
            float(actual),
            float(expected),
            rel_tol=relative_tolerance,
            abs_tol=absolute_tolerance,
        ):
            raise ReproductionError(
                f"numeric mismatch at {path}: {actual!r} != {expected!r}"
            )
        return
    if type(expected) is not type(actual):
        raise ReproductionError(
            f"type mismatch at {path}: {type(actual).__name__} != {type(expected).__name__}"
        )
    if isinstance(expected, dict):
        if set(expected) != set(actual):
            raise ReproductionError(f"mapping keys differ at {path}")
        for key in expected:
            compare_portable(
                expected[key],
                actual[key],
                path=f"{path}/{key}",
                relative_tolerance=relative_tolerance,
                absolute_tolerance=absolute_tolerance,
                compare_numbers=compare_numbers,
                numeric_overrides=numeric_overrides,
            )
        return
    if isinstance(expected, list):
        if len(expected) != len(actual):
            raise ReproductionError(f"list length differs at {path}")
        for index, (expected_item, actual_item) in enumerate(zip(expected, actual)):
            compare_portable(
                expected_item,
                actual_item,
                path=f"{path}/{index}",
                relative_tolerance=relative_tolerance,
                absolute_tolerance=absolute_tolerance,
                compare_numbers=compare_numbers,
                numeric_overrides=numeric_overrides,
            )
        return
    if expected != actual:
        raise ReproductionError(f"portable mismatch at {path}: {actual!r} != {expected!r}")


def resolve_pointer(value: Any, pointer: str) -> Any:
    if not pointer.startswith("/"):
        raise ReproductionError(f"invalid JSON pointer: {pointer}")
    current = value
    for raw_part in pointer.split("/")[1:]:
        part = raw_part.replace("~1", "/").replace("~0", "~")
        if isinstance(current, list):
            current = current[int(part)]
        elif isinstance(current, dict):
            current = current[part]
        else:
            raise ReproductionError(f"JSON pointer leaves a container: {pointer}")
    return current


def derived_numeric_overrides(
    expected: dict[str, Any], actual: dict[str, Any], policy: dict[str, Any]
) -> dict[str, tuple[float, float]]:
    """Propagate declared raw-error budgets through authenticated log-ratio identities.

    The budget is a release acceptance threshold, not a solver error estimate.
    Every error and order is recomputed independently within its own record;
    the asymmetric order interval follows monotonically from e +/- budget.
    """
    overrides: dict[str, tuple[float, float]] = {}

    def finite_number(value: Any, pointer: str) -> float:
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise ReproductionError(f"non-finite or non-numeric derived input at {pointer}")
        return float(value)

    def rounding_allowance(left: float, right: float, ulps: int) -> float:
        return ulps * max(math.ulp(left), math.ulp(right))

    for rule in policy.get("derived_quantities", []):
        if rule.get("kind") != "log2_absolute_error_ratio":
            raise ReproductionError("unknown derived-quantity comparison policy")
        budget = finite_number(rule["input_absolute_tolerance"], "/input_absolute_tolerance")
        identity_ulps = rule["identity_ulps"]
        interval_ulps = rule["interval_rounding_ulps"]
        if (budget <= 0 or type(identity_ulps) is not int or identity_ulps < 1
                or type(interval_ulps) is not int or interval_ulps < identity_ulps):
            raise ReproductionError("invalid derived-quantity error or rounding budget")
        for pointer in rule["family_pointers"]:
            records = []
            for label, value in (("expected", expected), ("actual", actual)):
                try:
                    family = resolve_pointer(value, pointer)
                    errors = family["absolute_errors"]
                    orders = family["observed_orders"]
                    lattice = family["lattice"]
                    continuum = finite_number(family["continuum"]["band_edge"], pointer + "/continuum/band_edge")
                except (KeyError, IndexError, TypeError) as exc:
                    raise ReproductionError(f"malformed derived family at {pointer}") from exc
                if (not isinstance(errors, list) or not isinstance(orders, list)
                        or not isinstance(lattice, list) or len(errors) < 2
                        or len(lattice) != len(errors) or len(orders) != len(errors) - 1):
                    raise ReproductionError(f"inconsistent derived family lengths at {pointer}")
                numeric_errors, numeric_orders = [], []
                for index, error in enumerate(errors):
                    error_pointer = f"{pointer}/absolute_errors/{index}"
                    error = finite_number(error, error_pointer)
                    if error <= budget:
                        raise ReproductionError(f"derived error is unresolved at the declared budget: {error_pointer}")
                    try:
                        gap = finite_number(lattice[index]["minimum_sampled_bloch_gap"], f"{pointer}/lattice/{index}/minimum_sampled_bloch_gap")
                    except (KeyError, TypeError) as exc:
                        raise ReproductionError(f"missing lattice gap at {pointer}") from exc
                    recomputed = abs(gap - continuum)
                    if abs(error - recomputed) > rounding_allowance(error, recomputed, identity_ulps):
                        raise ReproductionError(f"{label} error does not match its raw gap subtraction: {error_pointer}")
                    numeric_errors.append(error)
                for index, order in enumerate(orders):
                    order_pointer = f"{pointer}/observed_orders/{index}"
                    order = finite_number(order, order_pointer)
                    recomputed = math.log2(numeric_errors[index] / numeric_errors[index + 1])
                    if abs(order - recomputed) > rounding_allowance(order, recomputed, identity_ulps):
                        raise ReproductionError(f"{label} order does not match its error ratio: {order_pointer}")
                    numeric_orders.append(order)
                records.append((numeric_errors, numeric_orders))
            expected_errors, expected_orders = records[0]
            actual_errors, actual_orders = records[1]
            if len(expected_errors) != len(actual_errors):
                raise ReproductionError(f"derived family length changed at {pointer}")
            for index, (left, right) in enumerate(zip(expected_errors, actual_errors)):
                if abs(right - left) > budget:
                    raise ReproductionError(f"raw error exceeds its absolute portability budget: {pointer}/absolute_errors/{index}")
                overrides[f"{pointer}/absolute_errors/{index}"] = (0.0, budget)
            for index, order in enumerate(actual_orders):
                left, right = expected_errors[index:index + 2]
                lower = math.log2((left - budget) / (right + budget))
                upper = math.log2((left + budget) / (right - budget))
                cushion = interval_ulps * max(math.ulp(lower), math.ulp(upper),
                                             math.ulp(order), math.ulp(expected_orders[index]))
                if not lower - cushion <= order <= upper + cushion:
                    raise ReproductionError(f"order leaves its propagated log-ratio interval: {pointer}/observed_orders/{index}")
                allowed = max(abs(expected_orders[index] - lower),
                              abs(upper - expected_orders[index])) + cushion
                overrides[f"{pointer}/observed_orders/{index}"] = (0.0, allowed)
    return overrides


def validate_result(
    work: Path,
    expected_root: Path,
    step: dict[str, Any],
    actual_value: dict[str, Any],
    mode: str,
) -> bool:
    actual_path = work / step["output"]
    expected_path = expected_root / step["output"]
    actual_dynamic = validate_authenticated_inputs(work, actual_value, step)
    if mode == "exact":
        if actual_path.read_bytes() != expected_path.read_bytes():
            raise ReproductionError(f"exact byte mismatch: {step['output']}")
        return True

    expected_value = json.loads(expected_path.read_text(encoding="utf-8"))
    expected_dynamic = validate_authenticated_inputs(ROOT, expected_value, step, use_auxiliary=False)
    if actual_dynamic != expected_dynamic:
        raise ReproductionError("authenticated generated-result binding set changed")
    policy = step["comparison_policy"]
    policy_kind = str(policy.get("kind"))
    relative = float(policy.get("relative_tolerance", 0.0))
    absolute = float(policy.get("absolute_tolerance", 0.0))
    compare_portable(
        normalize_dynamic_hashes(expected_value, expected_dynamic),
        normalize_dynamic_hashes(actual_value, actual_dynamic),
        path="",
        relative_tolerance=relative,
        absolute_tolerance=absolute,
        compare_numbers=policy_kind in {"exact", "all_fields"},
        numeric_overrides=derived_numeric_overrides(expected_value, actual_value, policy),
    )
    for observable in step.get("headline_observables", []):
        pointer = str(observable["pointer"])
        compare_portable(
            resolve_pointer(expected_value, pointer),
            resolve_pointer(actual_value, pointer),
            path=pointer,
            relative_tolerance=float(observable["relative_tolerance"]),
            absolute_tolerance=float(observable["absolute_tolerance"]),
            compare_numbers=True,
        )
    return actual_path.read_bytes() == expected_path.read_bytes()


def execute(
    work: Path, manifest: dict[str, Any], *, mode: str, jobs: int
) -> dict[str, Any]:
    steps = {str(step["output"]): step for step in manifest["steps"]}
    pending = set(steps)
    complete: set[str] = set()
    running: dict[Future[tuple[dict[str, Any], float]], str] = {}
    durations: dict[str, float] = {}
    exact_matches = 0
    expected_root = work.parent / "expected"
    started = time.monotonic()

    with ThreadPoolExecutor(max_workers=jobs) as executor:
        while pending or running:
            ready = sorted(
                output
                for output in pending
                if set(steps[output]["dependencies"]) <= complete
            )
            while ready and len(running) < jobs:
                output = ready.pop(0)
                pending.remove(output)
                running[executor.submit(run_generator, work, steps[output])] = output
            if not running:
                raise ReproductionError("result graph stalled")
            finished, _ = wait(set(running), return_when=FIRST_COMPLETED)
            for future in finished:
                output = running.pop(future)
                value, elapsed = future.result()
                durations[output] = elapsed
                if validate_result(work, expected_root, steps[output], value, mode):
                    exact_matches += 1
                complete.add(output)
                print(
                    f"{len(complete):02d}/{len(steps)} {steps[output]['artifact_id']} "
                    f"{elapsed:.2f}s",
                    flush=True,
                )

    frontier = json.loads((work / manifest["frontier_output"]).read_text(encoding="utf-8"))
    return {
        "schema": "NSC-PUBLIC-REPRODUCTION-SUMMARY-v1",
        "release_version": manifest["release_version"],
        "release_spec_sha256": manifest["release_spec_sha256"],
        "release_source_commit": manifest["release_source_commit"],
        "historical_result_count": manifest["historical_result_count"],
        "scoped_result_count": len(manifest["scoped_follow_up_outputs"]),
        "mode": mode,
        "jobs": jobs,
        "result_count": len(complete),
        "exact_byte_matches": exact_matches,
        "portable_matches": len(complete),
        "elapsed_seconds": time.monotonic() - started,
        "frontier_artifact_id": frontier["artifact_id"],
        "frontier_classification": frontier["classification"],
        "slowest_steps": [
            {"output": output, "seconds": seconds}
            for output, seconds in sorted(
                durations.items(), key=lambda item: item[1], reverse=True
            )[:5]
        ],
    }


def main() -> int:
    args = parse_args()
    jobs = resolve_jobs(args.jobs)
    manifest = load_manifest()
    validate_checkout(manifest)
    workspace = Path(tempfile.mkdtemp(prefix="nsc-public-reproduction-"))
    success = False
    try:
        work = prepare_workspace(workspace, manifest)
        summary = execute(work, manifest, mode=args.mode, jobs=jobs)
        summary["workspace_retained"] = bool(args.keep_workspace)
        if args.keep_workspace:
            summary["workspace"] = str(workspace)
        rendered = json.dumps(summary, indent=2, sort_keys=True) + "\n"
        print(rendered, end="")
        if args.json_summary is not None:
            args.json_summary.parent.mkdir(parents=True, exist_ok=True)
            args.json_summary.write_text(rendered, encoding="utf-8")
        success = True
        return 0
    except Exception as exc:
        print(f"reproduction failed; workspace retained at {workspace}", file=sys.stderr)
        print(str(exc), file=sys.stderr)
        return 1
    finally:
        if success and not args.keep_workspace:
            shutil.rmtree(workspace)


if __name__ == "__main__":
    raise SystemExit(main())
