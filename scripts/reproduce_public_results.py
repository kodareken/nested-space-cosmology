#!/usr/bin/env python3
"""Recompute the curated result graph without modifying the checkout."""

from __future__ import annotations

import argparse
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
    if not isinstance(steps, list) or len(steps) != 59:
        raise ReproductionError("public result manifest must contain 59 steps")
    if historical != 58:
        raise ReproductionError("historical result count must remain 58")
    return value


def validate_checkout(manifest: dict[str, Any]) -> None:
    outputs: set[str] = set()
    for step in manifest["steps"]:
        output = str(step["output"])
        generator = str(step["generator"])
        if output in outputs:
            raise ReproductionError(f"duplicate result output: {output}")
        outputs.add(output)
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
    process = subprocess.run(
        [sys.executable, str(work / step["generator"])],
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
    if raw != canonical_json(value):
        raise ReproductionError(f"result is not canonical JSON: {step['output']}")
    if value.get("artifact_id") != step["artifact_id"] or value.get("terminal") is not True:
        raise ReproductionError(f"invalid result identity or terminal: {step['output']}")
    return value, elapsed


def authenticated_entries(value: dict[str, Any]):
    for key, collection in value.items():
        if not key.startswith("authenticated"):
            continue
        if isinstance(collection, dict) and {"path", "sha256"} <= set(collection):
            yield collection
            continue
        if isinstance(collection, dict):
            entries = collection.values()
        elif isinstance(collection, list):
            entries = collection
        else:
            raise ReproductionError(f"malformed {key} collection")
        for entry in entries:
            if not isinstance(entry, dict) or not {"path", "sha256"} <= set(entry):
                raise ReproductionError(f"malformed entry in {key}")
            yield entry


def auxiliary_lookup(work: Path, step: dict[str, Any] | None) -> dict[str, Path]:
    if not step:
        return {}
    root = work / "provenance" / Path(str(step["output"])).stem
    lookup: dict[str, Path] = {}
    for auxiliary in step.get("auxiliary_inputs", []):
        relative = str(auxiliary["path"])
        lookup[relative] = root / relative
    return lookup


def validate_authenticated_inputs(
    work: Path, value: dict[str, Any], step: dict[str, Any] | None = None
) -> None:
    lookup = auxiliary_lookup(work, step)
    for entry in authenticated_entries(value):
        relative = entry.get("path")
        expected = entry.get("sha256")
        if not isinstance(relative, str) or not isinstance(expected, str):
            raise ReproductionError("malformed authenticated input entry")
        path = lookup.get(relative, work / relative)
        if not path.is_file() or sha256(path) != expected:
            raise ReproductionError(f"authenticated input mismatch: {relative}")
        artifact_id = entry.get("artifact_id")
        if artifact_id is not None and relative.startswith("results/"):
            actual_id = json.loads(path.read_text(encoding="utf-8")).get("artifact_id")
            if actual_id != artifact_id:
                raise ReproductionError(f"authenticated artifact mismatch: {relative}")


def normalize_dynamic_hashes(value: dict[str, Any]) -> dict[str, Any]:
    normalized = json.loads(json.dumps(value))
    for entry in authenticated_entries(normalized):
        path = entry.get("path", "")
        if isinstance(path, str) and path.startswith("results/"):
            entry["sha256"] = "<validated-generated-result>"
    return normalized


def compare_portable(
    expected: Any,
    actual: Any,
    *,
    path: str,
    relative_tolerance: float,
    absolute_tolerance: float,
    compare_numbers: bool,
) -> None:
    if isinstance(expected, bool) or isinstance(actual, bool):
        if expected is not actual:
            raise ReproductionError(f"portable mismatch at {path}: {actual!r} != {expected!r}")
        return
    if isinstance(expected, (int, float)) and isinstance(actual, (int, float)):
        if not compare_numbers:
            return
        if isinstance(expected, int) and isinstance(actual, int):
            if expected != actual:
                raise ReproductionError(f"integer mismatch at {path}: {actual} != {expected}")
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


def validate_result(
    work: Path,
    expected_root: Path,
    step: dict[str, Any],
    actual_value: dict[str, Any],
    mode: str,
) -> bool:
    actual_path = work / step["output"]
    expected_path = expected_root / step["output"]
    validate_authenticated_inputs(work, actual_value, step)
    if mode == "exact":
        if actual_path.read_bytes() != expected_path.read_bytes():
            raise ReproductionError(f"exact byte mismatch: {step['output']}")
        return True

    expected_value = json.loads(expected_path.read_text(encoding="utf-8"))
    policy = step["comparison_policy"]
    policy_kind = str(policy.get("kind"))
    relative = float(policy.get("relative_tolerance", 0.0))
    absolute = float(policy.get("absolute_tolerance", 0.0))
    compare_portable(
        normalize_dynamic_hashes(expected_value),
        normalize_dynamic_hashes(actual_value),
        path="",
        relative_tolerance=relative,
        absolute_tolerance=absolute,
        compare_numbers=policy_kind in {"exact", "all_fields"},
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
