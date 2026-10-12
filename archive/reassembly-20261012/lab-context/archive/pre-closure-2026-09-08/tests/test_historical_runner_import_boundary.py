"""Prevent new dependencies on the frozen versioned GR-0 runner chain."""

from __future__ import annotations

import ast
from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]
VERSIONED = re.compile(r"run_fgc_gr0_calibration_v[0-9]+\Z")
SCAN_ROOTS = (ROOT / "src", ROOT / "scripts", ROOT / "tests")

# These files are immutable historical consumers or tests of that exact chain.
# New files must depend on stable runtime modules instead.
GRANDFATHERED = frozenset(
    {
        "scripts/reproduce_fgc_cal3_pref5.py",
        "scripts/reproduce_fgc_cal4_pref6.py",
        "scripts/reproduce_fgc_hlt5_mon5.py",
        "scripts/run_fgc_gr0_calibration_v10.py",
        "scripts/run_fgc_gr0_calibration_v11.py",
        "scripts/run_fgc_gr0_calibration_v12.py",
        "scripts/run_fgc_gr0_calibration_v13.py",
        "scripts/run_fgc_gr0_calibration_v14.py",
        "scripts/run_fgc_gr0_calibration_v7.py",
        "scripts/run_fgc_gr0_calibration_v8.py",
        "scripts/run_fgc_gr0_calibration_v9.py",
        "scripts/run_fgc_rsp1_resolution_study.py",
        "scripts/run_fgc_rsp2_constraint_study.py",
        "scripts/run_fgc_src2_affine_capture.py",
        "src/recursive_horizons/fgc/evolution/proto19_progression_inputs.py",
        "tests/test_fgc_gr0_campaign_runner_v10.py",
        "tests/test_fgc_gr0_campaign_runner_v11.py",
        "tests/test_fgc_gr0_campaign_runner_v12.py",
        "tests/test_fgc_gr0_campaign_runner_v13.py",
        "tests/test_fgc_gr0_campaign_runner_v14.py",
        "tests/test_fgc_gr0_campaign_runner_v6.py",
        "tests/test_fgc_gr0_campaign_runner_v7.py",
        "tests/test_fgc_gr0_campaign_runner_v8.py",
        "tests/test_fgc_gr0_campaign_runner_v9.py",
        "tests/test_fgc_proto19_gr0_static_factory.py",
        "tests/test_fgc_proto19_progression_inputs.py",
    }
)


def _versioned_imports(path: Path) -> tuple[str, ...]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    imports: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module == "scripts":
            imports.extend(
                alias.name
                for alias in node.names
                if VERSIONED.fullmatch(alias.name)
            )
        elif (
            isinstance(node, ast.ImportFrom)
            and node.module is not None
            and node.module.startswith("scripts.")
        ):
            name = node.module.rsplit(".", 1)[-1]
            if VERSIONED.fullmatch(name):
                imports.append(name)
        elif isinstance(node, ast.Import):
            for alias in node.names:
                name = alias.name.rsplit(".", 1)[-1]
                if alias.name.startswith("scripts.") and VERSIONED.fullmatch(name):
                    imports.append(name)
    return tuple(sorted(imports))


class HistoricalRunnerImportBoundaryTests(unittest.TestCase):
    def test_only_grandfathered_files_import_versioned_runner_modules(self) -> None:
        observed: dict[str, tuple[str, ...]] = {}
        for directory in SCAN_ROOTS:
            for path in sorted(directory.rglob("*.py")):
                relative = path.relative_to(ROOT).as_posix()
                imports = _versioned_imports(path)
                if imports:
                    observed[relative] = imports
        self.assertEqual(set(observed), set(GRANDFATHERED))
        self.assertTrue(all(observed.values()))

    def test_grandfathered_paths_exist_and_are_outside_active_extension_names(self) -> None:
        for relative in GRANDFATHERED:
            self.assertTrue((ROOT / relative).is_file(), relative)
        self.assertFalse(
            any(
                any(token in relative for token in ("sgb1", "def1", "col1", "rob1", "run2"))
                for relative in GRANDFATHERED
            )
        )


if __name__ == "__main__":
    unittest.main()
