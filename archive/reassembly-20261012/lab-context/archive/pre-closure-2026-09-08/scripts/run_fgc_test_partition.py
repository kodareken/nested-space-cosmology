#!/usr/bin/env python3
"""Run the repository tests under their hash-owning arithmetic backend."""

from __future__ import annotations

import argparse
import platform
from pathlib import Path
import subprocess
import sys
import tomllib
import unittest
from typing import Any, Iterable, Iterator

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
CONFIG = REPOSITORY / "configs/fgc/fgc-runtime-matrix.toml"
TESTS = REPOSITORY / "tests"
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

# These unmodified tests require a finite decimal-conversion cap. Historical
# exact-certificate imports intentionally disable that process-global cap.
# Keep both owners intact and run the finite-cap domain in a fresh interpreter.
PROCESS_ISOLATED_MODULES = {
    "tests.test_fgc_tdg11_imp1_enclosure": 4300,
}
PROCESS_ISOLATION_TIMEOUT_SECONDS = 600


def _load_config() -> dict[str, Any]:
    with CONFIG.open("rb") as handle:
        value = tomllib.load(handle)
    if set(value) != {
        "schema_version", "artifact_id", "project_version", "analytic",
        "evolution_restart", "partition", "paper",
    }:
        raise ValueError("FGC runtime-matrix root differs")
    if (
        value["schema_version"] != 1
        or value["artifact_id"] != "FGC-RUNTIME-MATRIX-1"
        or value["project_version"] != "0.11.0"
    ):
        raise ValueError("FGC runtime-matrix identity differs")
    partition = value["partition"]
    if not isinstance(partition, dict) or set(partition) != {
        "runtime_bound_test_modules",
        "raw_bundle_paths",
        "complete_raw_absence_skips_only_runtime_bound_tests",
        "partial_raw_presence_fails_closed",
        "proto14_postlaunch",
        "proto18_pref27_postlaunch",
        "tdg9_ur1_pref1_postlaunch",
    }:
        raise ValueError("FGC runtime-matrix partition differs")
    postlaunch = partition.get("proto14_postlaunch")
    if not isinstance(postlaunch, dict) or set(postlaunch) != {
        "cal11_pref22_reproducer", "cal11_pref22_result",
        "sealed_prelaunch_test_ids",
    }:
        raise ValueError("FGC PROTO14 postlaunch routing differs")
    if (
        postlaunch["cal11_pref22_reproducer"]
        != "scripts/reproduce_fgc_cal11_pref22.py"
        or postlaunch["cal11_pref22_result"]
        != "results/fgc-1-cal11-pref22.json"
        or not isinstance(postlaunch["sealed_prelaunch_test_ids"], list)
        or len(postlaunch["sealed_prelaunch_test_ids"]) != 3
        or len(set(postlaunch["sealed_prelaunch_test_ids"])) != 3
        or any(
            not isinstance(test_id, str)
            for test_id in postlaunch["sealed_prelaunch_test_ids"]
        )
    ):
        raise ValueError("FGC PROTO14 postlaunch routing contract differs")
    pref27_postlaunch = partition.get("proto18_pref27_postlaunch")
    if (
        not isinstance(pref27_postlaunch, dict)
        or set(pref27_postlaunch) != {"sealed_prelaunch_test_ids"}
    ):
        raise ValueError("FGC PREF27 postlaunch routing differs")
    pref27_sealed = pref27_postlaunch["sealed_prelaunch_test_ids"]
    if (
        not isinstance(pref27_sealed, list)
        or len(pref27_sealed) != 13
        or len(set(pref27_sealed)) != 13
        or any(not isinstance(test_id, str) for test_id in pref27_sealed)
    ):
        raise ValueError("FGC PREF27 postlaunch routing contract differs")
    ur1_postlaunch = partition.get("tdg9_ur1_pref1_postlaunch")
    if not isinstance(ur1_postlaunch, dict) or set(ur1_postlaunch) != {
        "ur1_pref1_reproducer",
        "sealed_prelaunch_test_ids",
        "evolution_restart_reconstruction_test_ids",
    }:
        raise ValueError("FGC UR1-PREF1 postlaunch routing differs")
    ur1_sealed = ur1_postlaunch["sealed_prelaunch_test_ids"]
    ur1_routed = ur1_postlaunch["evolution_restart_reconstruction_test_ids"]
    if (
        ur1_postlaunch["ur1_pref1_reproducer"]
        != "scripts/reproduce_fgc_tdg9_ur1_pref1.py"
        or not isinstance(ur1_sealed, list)
        or len(ur1_sealed) != 9
        or len(set(ur1_sealed)) != 9
        or any(
            not isinstance(test_id, str)
            or not test_id.startswith("tests.test_")
            or test_id != test_id.strip()
            or not test_id.split(".")[-1].startswith("test_")
            for test_id in ur1_sealed
        )
        or not isinstance(ur1_routed, list)
        or len(ur1_routed) != 2
        or len(set(ur1_routed)) != 2
        or any(
            not isinstance(test_id, str)
            or not test_id.startswith("tests.test_")
            or test_id != test_id.strip()
            or not test_id.split(".")[-1].startswith("test_")
            for test_id in ur1_routed
        )
        or set(ur1_sealed) & set(ur1_routed)
    ):
        raise ValueError("FGC UR1-PREF1 postlaunch routing contract differs")
    bound = set(partition["runtime_bound_test_modules"])
    owned_ids = list(ur1_sealed) + list(ur1_routed)
    if any(test_id.split(".")[1] in bound for test_id in owned_ids):
        raise ValueError(
            "FGC UR1-PREF1 routed tests are already runtime-bound at module grain"
        )
    combined_ids = (
        list(postlaunch["sealed_prelaunch_test_ids"])
        + list(pref27_sealed)
        + owned_ids
    )
    if len(combined_ids) != len(set(combined_ids)):
        raise ValueError("FGC postlaunch test IDs are not unique across routing sections")
    return value


def _environment() -> dict[str, str]:
    configuration = np.show_config(mode="dicts")
    dependencies = configuration.get("Build Dependencies", {})
    blas = dependencies.get("blas", {}) if isinstance(dependencies, dict) else {}
    if not isinstance(blas, dict):
        raise ValueError("NumPy BLAS configuration is unavailable")
    return {
        "python_version": platform.python_version(),
        "python_implementation": platform.python_implementation(),
        "numpy_version": np.__version__,
        "blas_name": str(blas.get("name")),
        "blas_version": str(blas.get("version")),
        "system": platform.system(),
        "machine": platform.machine(),
    }


def _validate_environment(config: dict[str, Any], partition: str) -> None:
    expected = config["analytic" if partition == "analytic" else "evolution_restart"]
    observed = _environment()
    if observed != expected:
        raise ValueError(
            f"{partition} test runtime differs: expected {expected!r}, "
            f"observed {observed!r}, executable={sys.executable}"
        )


def _module_names() -> tuple[str, ...]:
    return tuple(
        f"tests.{path.stem}"
        for path in sorted(TESTS.glob("test*.py"))
        if path.is_file()
    )


def _select_modules(
    modules: Iterable[str], runtime_bound: set[str], partition: str
) -> tuple[str, ...]:
    if partition == "analytic":
        return tuple(
            module for module in modules if module.rsplit(".", 1)[-1] not in runtime_bound
        )
    return tuple(
        module for module in modules if module.rsplit(".", 1)[-1] in runtime_bound
    )


def _raw_bundle_state(config: dict[str, Any]) -> str:
    paths = [REPOSITORY / item for item in config["partition"]["raw_bundle_paths"]]
    present = [path.is_file() for path in paths]
    if all(present):
        return "complete"
    if not any(present):
        return "absent"
    missing = [path.relative_to(REPOSITORY).as_posix() for path, ok in zip(paths, present) if not ok]
    raise ValueError(f"FGC runtime-bound raw bundle is partial; missing {missing}")


def _iter_cases(suite: unittest.TestSuite) -> Iterator[unittest.TestCase]:
    """Flatten one loaded suite without changing its test identities."""

    for item in suite:
        if isinstance(item, unittest.TestSuite):
            yield from _iter_cases(item)
        else:
            yield item


def _split_process_isolated_cases(
    suite: unittest.TestSuite,
) -> tuple[unittest.TestSuite, dict[str, tuple[str, ...]]]:
    """Route exact existing test IDs without skipping or replacing a test."""

    retained = []
    isolated: dict[str, list[str]] = {}
    cases = tuple(_iter_cases(suite))
    for case in cases:
        module = type(case).__module__
        if module in PROCESS_ISOLATED_MODULES:
            test_id = case.id()
            if not test_id.startswith(module + "."):
                raise ValueError("isolated test identity differs from its module")
            isolated.setdefault(module, []).append(test_id)
        else:
            retained.append(case)
    if len(retained) + sum(map(len, isolated.values())) != len(cases):
        raise RuntimeError("process isolation lost a selected test")
    return unittest.TestSuite(retained), {
        module: tuple(test_ids) for module, test_ids in isolated.items()
    }


def _run_process_isolated_cases(groups: dict[str, tuple[str, ...]]) -> int:
    """Run each process-state domain under the same Python/NumPy installation."""

    child = (
        "from pathlib import Path; import sys, unittest; "
        "root=Path(sys.argv[1]); "
        "sys.path[:0]=[str(root),str(root/'src')]; "
        "names=sys.argv[2:]; "
        "suite=unittest.defaultTestLoader.loadTestsFromNames(names); "
        "assert suite.countTestCases()==len(names), 'isolated case count differs'; "
        "result=unittest.TextTestRunner(verbosity=2).run(suite); "
        "raise SystemExit(0 if result.wasSuccessful() else 1)"
    )
    count = 0
    for module, test_ids in groups.items():
        if (
            module not in PROCESS_ISOLATED_MODULES
            or type(test_ids) is not tuple
            or not test_ids
            or len(set(test_ids)) != len(test_ids)
            or any(type(name) is not str or not name.startswith(module + ".")
                   for name in test_ids)
        ):
            raise ValueError("isolated process group differs from its declared owner")
        limit = PROCESS_ISOLATED_MODULES[module]
        print(
            f"FGC process isolation module={module} cases={len(test_ids)} "
            f"int_max_str_digits={limit}",
            flush=True,
        )
        completed = subprocess.run(
            [sys.executable, "-I", "-B", "-X", f"int_max_str_digits={limit}",
             "-c", child, str(REPOSITORY), *test_ids],
            cwd=REPOSITORY,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=PROCESS_ISOLATION_TIMEOUT_SECONDS,
            check=False,
        )
        if completed.stdout:
            print(completed.stdout, end="", flush=True)
        if completed.stderr:
            print(completed.stderr, end="", file=sys.stderr, flush=True)
        if completed.returncode != 0:
            raise RuntimeError(f"isolated process tests failed: {module}")
        count += len(test_ids)
    return count


def _run_ur1_pref1_check(config: dict[str, Any]) -> None:
    """Bind compact UR1-PREF1 evidence before postlaunch analytic routing."""

    postlaunch = config["partition"]["tdg9_ur1_pref1_postlaunch"]
    reproducer = REPOSITORY / postlaunch["ur1_pref1_reproducer"]
    if not reproducer.is_file():
        raise ValueError("UR1-PREF1 postlaunch verifier is absent")
    completed = subprocess.run(
        [sys.executable, str(reproducer)],
        cwd=REPOSITORY,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if completed.returncode != 0:
        detail = completed.stderr.strip() or completed.stdout.strip()
        raise RuntimeError(f"UR1-PREF1 check-only verification failed: {detail}")
    print(
        "UR1-PREF1 check-only/store-blind verification passed before postlaunch sealing",
        flush=True,
    )


def _run_cal11_pref22_check(config: dict[str, Any]) -> None:
    """Bind complete PROTO14 raw evidence before postlaunch test routing."""

    postlaunch = config["partition"]["proto14_postlaunch"]
    reproducer = REPOSITORY / postlaunch["cal11_pref22_reproducer"]
    result = REPOSITORY / postlaunch["cal11_pref22_result"]
    if not reproducer.is_file() or not result.is_file():
        raise ValueError("CAL11-PREF22 postlaunch verifier is absent")
    completed = subprocess.run(
        [sys.executable, str(reproducer), "--check", "--output", str(result)],
        cwd=REPOSITORY,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if completed.returncode != 0:
        detail = completed.stderr.strip() or completed.stdout.strip()
        raise RuntimeError(f"CAL11-PREF22 check-only verification failed: {detail}")
    print("CAL11-PREF22 check-only verification passed before postlaunch sealing", flush=True)


def _seal_prelaunch_tests(
    suite: unittest.TestSuite,
    sealed_ids: Iterable[str],
    *,
    seal_id: str = "FGC-1-CAL11-PREF22",
) -> unittest.TestSuite:
    """Remove exactly the sealed prelaunch tests from a loaded partition suite."""

    sealed = tuple(sealed_ids)
    cases = tuple(_iter_cases(suite))
    available = {case.id() for case in cases}
    missing = set(sealed) - available
    if missing:
        raise ValueError(
            "configured postlaunch prelaunch-only test IDs are absent: "
            f"{sorted(missing)}"
        )
    retained = tuple(case for case in cases if case.id() not in sealed)
    filtered = {case.id() for case in cases} - {case.id() for case in retained}
    if filtered != set(sealed) or len(cases) - len(retained) != len(sealed):
        raise RuntimeError("postlaunch routing filtered tests outside the sealed set")
    print(
        f"FGC postlaunch seal={seal_id} cases={len(sealed)}",
        flush=True,
    )
    for test_id in sealed:
        print(f"SEALED prelaunch-only test: {test_id}", flush=True)
    return unittest.TestSuite(retained)


def _analytic_module_name(test_id: str) -> str:
    parts = test_id.split(".")
    if len(parts) < 4 or parts[0] != "tests":
        raise ValueError(f"reconstruction test ID is invalid: {test_id!r}")
    return ".".join(parts[:2])


def _route_reconstruction_cases(
    source_suite: unittest.TestSuite,
    routed_ids: Iterable[str],
    *,
    owner: str = "FGC-1-TDG9-UR1-PREF1",
) -> tuple[unittest.TestSuite, tuple[unittest.TestCase, ...]]:
    """Prove reconstruction IDs exist in analytic and extract exactly those cases."""

    routed = tuple(routed_ids)
    cases = tuple(_iter_cases(source_suite))
    available = {case.id(): case for case in cases}
    missing = [test_id for test_id in routed if test_id not in available]
    if missing:
        raise ValueError(
            "configured reconstruction test IDs are absent from the analytic suite: "
            f"{missing}"
        )
    extracted = tuple(available[test_id] for test_id in routed)
    retained = tuple(case for case in cases if case.id() not in set(routed))
    filtered = {case.id() for case in cases} - {case.id() for case in retained}
    if (
        filtered != set(routed)
        or len(extracted) != len(routed)
        or len(cases) - len(retained) != len(routed)
        or tuple(case.id() for case in extracted) != routed
    ):
        raise RuntimeError(
            "reconstruction routing filtered tests outside the routed set"
        )
    print(
        f"FGC reconstruction route owner={owner} cases={len(routed)}",
        flush=True,
    )
    for test_id in routed:
        print(f"ROUTED reconstruction test: {test_id}", flush=True)
    return unittest.TestSuite(retained), extracted


def _append_reconstruction_cases(
    suite: unittest.TestSuite,
    cases: Iterable[unittest.TestCase],
) -> unittest.TestSuite:
    """Append exact analytic reconstruction cases to the evolution suite."""

    incoming = tuple(cases)
    added_ids = tuple(case.id() for case in incoming)
    if len(set(added_ids)) != len(added_ids):
        raise ValueError("reconstruction cases are not unique")
    retained = tuple(_iter_cases(suite))
    existing_ids = {case.id() for case in retained}
    overlap = existing_ids.intersection(added_ids)
    if overlap:
        raise ValueError(
            "reconstruction cases already present in evolution-restart: "
            f"{sorted(overlap)}"
        )
    combined = unittest.TestSuite(retained + incoming)
    combined_ids = {case.id() for case in _iter_cases(combined)}
    if (
        combined.countTestCases() != len(retained) + len(incoming)
        or combined_ids - existing_ids != set(added_ids)
    ):
        raise RuntimeError("reconstruction append changed unlisted tests")
    print(
        f"FGC reconstruction append partition=evolution-restart cases={len(incoming)}",
        flush=True,
    )
    return combined


def _analytic_reconstruction_cases(
    config: dict[str, Any],
    runtime_bound: set[str],
) -> tuple[unittest.TestCase, ...]:
    """Load exactly the two analytic reconstruction cases for evolution-restart."""

    routed = tuple(
        config["partition"]["tdg9_ur1_pref1_postlaunch"][
            "evolution_restart_reconstruction_test_ids"
        ]
    )
    analytic_modules = set(
        _select_modules(_module_names(), runtime_bound, "analytic")
    )
    source_modules: list[str] = []
    for test_id in routed:
        module = _analytic_module_name(test_id)
        if module not in analytic_modules:
            raise ValueError(
                f"reconstruction test is not in the analytic suite: {test_id}"
            )
        if module not in source_modules:
            source_modules.append(module)
    source = unittest.defaultTestLoader.loadTestsFromNames(tuple(source_modules))
    _retained, extracted = _route_reconstruction_cases(source, routed)
    if len(extracted) != 2:
        raise RuntimeError(
            "evolution reconstruction routing did not return exactly two cases"
        )
    return extracted


def _apply_tdg9_ur1_pref1_analytic_routing(
    config: dict[str, Any],
    suite: unittest.TestSuite,
) -> tuple[unittest.TestSuite, tuple[unittest.TestCase, ...]]:
    """Store-blind compact check, then seal nine tests and extract two reconstructions."""

    _run_ur1_pref1_check(config)
    postlaunch = config["partition"]["tdg9_ur1_pref1_postlaunch"]
    sealed = _seal_prelaunch_tests(
        suite,
        postlaunch["sealed_prelaunch_test_ids"],
        seal_id="FGC-1-TDG9-UR1-PREF1",
    )
    return _route_reconstruction_cases(
        sealed,
        postlaunch["evolution_restart_reconstruction_test_ids"],
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--partition", choices=("analytic", "evolution-restart"), required=True
    )
    parser.add_argument(
        "--portable-clean-clone",
        action="store_true",
        help=(
            "run portable clean-clone validation without claiming canonical "
            "local arithmetic reproduction"
        ),
    )
    args = parser.parse_args()
    config = _load_config()
    runtime_bound = set(config["partition"]["runtime_bound_test_modules"])
    if args.partition == "evolution-restart":
        raw_state = _raw_bundle_state(config)
        if raw_state == "absent":
            print(
                "SKIP evolution-restart test partition: complete ignored raw bundle is "
                "absent; the repository audit will verify compact hash-bound evidence"
            )
            return
        if args.portable_clean_clone:
            raise ValueError(
                "portable clean-clone mode cannot consume the local runtime-bound raw bundle"
            )
    if args.portable_clean_clone:
        print(
            "PORTABLE clean-clone validation: runtime fingerprint is reported but is not "
            "canonical arithmetic provenance",
            flush=True,
        )
    else:
        _validate_environment(config, args.partition)
    modules = _select_modules(_module_names(), runtime_bound, args.partition)
    if not modules:
        raise ValueError(f"{args.partition} test partition is empty")
    suite = unittest.defaultTestLoader.loadTestsFromNames(modules)
    if args.partition == "analytic":
        suite = _seal_prelaunch_tests(
            suite,
            config["partition"]["proto18_pref27_postlaunch"][
                "sealed_prelaunch_test_ids"
            ],
            seal_id="FGC-1-PRO18-PREF27",
        )
        suite, _routed = _apply_tdg9_ur1_pref1_analytic_routing(config, suite)
    if args.partition == "evolution-restart" and raw_state == "complete":
        _run_cal11_pref22_check(config)
        suite = _seal_prelaunch_tests(
            suite,
            config["partition"]["proto14_postlaunch"]["sealed_prelaunch_test_ids"],
        )
        suite = _append_reconstruction_cases(
            suite,
            _analytic_reconstruction_cases(config, runtime_bound),
        )
    count = suite.countTestCases()
    suite, isolated = _split_process_isolated_cases(suite)
    isolated_count = sum(map(len, isolated.values()))
    if suite.countTestCases() + isolated_count != count:
        raise RuntimeError("process-isolated partition count differs")
    print(
        f"FGC test partition={args.partition} cases={count} "
        f"in_process={suite.countTestCases()} isolated={isolated_count} "
        f"python={platform.python_version()} numpy={np.__version__}",
        flush=True,
    )
    if _run_process_isolated_cases(isolated) != isolated_count:
        raise RuntimeError("process-isolated execution count differs")
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if not result.wasSuccessful():
        raise SystemExit(1)


if __name__ == "__main__":
    main()
